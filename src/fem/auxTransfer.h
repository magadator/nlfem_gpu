#include <vector>
using namespace std;
//struct geomF_2_femF {
//  vector<int> closestFEMnode;
//};

struct auxStencil {
  //normal direction from FEM
  double sign;
  //displacement transfer
  vector<int>    nodes;
  vector<double> coefs;
  vector<vector<int> > closestGeomNode;
  //load transfer
  vector<int> closestFEMnode;
  vector<vector<double> > loadCoefs;
};

