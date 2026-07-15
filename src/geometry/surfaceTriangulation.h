#ifndef SURFTRI
#define SURFTRI
//
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
using namespace std;
namespace KARMA { 

class surfaceTriangulation {
public:
void smartTriangulationRead(vector<double>& nodes,                //xyz
			      vector<int>&    faces,                //connectivity
			      string&         filename,             //VTK filename
			      int&            globalNodeNum,        //global number of unqiue nodes
			      int&            globalFaceNum,        //global number of elements
			      vector<int>&    elementTypes,         //which element types were specified in input file
			      vector<int>&    elementTypeVector,    //tag each global element ID with its element type
			      vector<int>&    ndofPerNode,          //unique list of ndof per node
			      vector<int>&    componentIDs,         //element component IDs
			      bool&           revertNormalsByComp); //should we revert normals on the listed comp IDs

void readTriangulation(vector<double>& nodes,
			 vector<int>&    faces,
			 string&         filename,
			 int&            globalNodeNum,
			 int&            globalFaceNum,
			 int&            nodePerElement);
void writeTriangulation_parallel_unformatted(string& filename,
					     int& globalNodeCount,
					     int& globalFaceCount,
					     vector<double>& localNodeForces,
					     vector<double>& localNodeDispls,
					     map<int,int>& nodeLocal2Global,
					     map<int,int>& elemLocal2Global,
					     int& localNodeCount,
					     int& localFaceCount,
					     vector<double>& nodes,
					     vector<int>& faces,
					     vector<int>& nnodesPerElement,
					     vector<double>& vonMises);
    void readTriangulation_AUX(vector<double>& nodes,
                               vector<int>&    faces,
                               string&         filename,
                               int&            globalNodeNum,
                               int&            globalFaceNum,
                               int&            nodePerElement,
                               vector<int>&    isauxnode);
    void readTriangulation_AUX_FSI(vector<double>& nodes,
                                   string&         filename);
    void writeTriangulation(string filename,
                            vector<double>& Nodes,
                            vector<int>&    Faces,
                            int& DIM,
                            int& nodeNum,
                            int &faceNum);    
    void smartWriteTriangulation(string filename,
                                 vector<double>& Nodes,
                                 vector<int>&    Faces,
                                 int& DIM,
                                 int& nodeNum,
                                 int &faceNum,
                                 int& numberOfElementTypes,
                                 vector<int>& nnodesPerElement);
    void smartWriteTriangulation_withForces(string filename,
                                            vector<double>& Nodes,
                                            vector<int>&    Faces,
                                            int& DIM,
                                            int& nodeNum,
                                            int &faceNum,
                                            int& numberOfElementTypes,
                                            vector<int>& nnodesPerElement,
                                            vector<double>& forces,
                                            vector<int>& orderedList,
					    vector<double>& disps,
					    vector<double>& vonMises);
    void read_preStrained_nodeLocs(vector<vector<double> >& nodes,
                                   string filename);
    void writeTriangulation_data(string filename,
                                 vector<vector<double> >& Nodes,
                                 vector<int>& Faces,
                                 int& DIM,
                                 vector<vector<double> >& forces,
                                 vector<vector<double> >& disp,
                                 vector<vector<double> >& vel,
                                 vector<double>& data,
                                 vector<vector<double> >& Vnt,
                                 vector<vector<double> >& V1t,
                                 vector<vector<double> >& V2t);
    void writeTriangulation_aux(string filename,
                                vector<double>& Nodes,
                                vector<int>& Faces,
                                int& DIM);                             
    void writeVTK(string filename,
                  vector<vector<double> >& Nodes,
                  vector<int>& Faces,
                  vector<vector<double> >& Data,
                  vector<string>& varnames);
    void writeVTK_T(string filename,
                    vector<vector<double> >& Nodes,
                    vector<int>& Faces,
                    vector<vector<double> >& Data,
                    vector<string>& varnames);
    void writeVTK_disp(string filename,
                       vector<vector<double> >& Nodes,
                       vector<int>& Faces,
                       vector<vector<double> >& disps);
    void writeVTK_withVectors(string filename,
                              vector<vector<double> >& Nodes,
                              vector<int>& Faces,
                              vector<vector<double> >& Vn,
                              vector<vector<double> >& V1,
                              vector<vector<double> >& V2);
    vector<vector<double> > nodes_;
    vector<vector<int> >    faces_;
    //temp aux geometry vars
    vector<vector<double> > auxNodes_;
    vector<vector<int> >    auxFaces_;
    int auxNodeNum_;
    int auxFaceNum_;    
    //    
    vector<vector<double> > Vnt_surf_;
    vector<int>        auxComponents_;
    vector<int>        components_;
    string surfaceFilename_;
    int globalNodeNum_;
    int globalFaceNum_; 
    //
  };
  void writeTriangle(string filename,vector<double>& Node1,vector<double>& Node2,vector<double>& Node3);
}
#endif
