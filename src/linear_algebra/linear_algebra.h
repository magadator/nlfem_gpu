#ifndef KARMA_LINEAR_ALGEBRA_H
#define KARMA_LINEAR_ALGEBRA_H
//
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include "mpi.h"

using namespace std;

namespace KARMA {
  //
  double checksign(double& num);
  //
  void matVec_mult(vector<vector<double> > &a,
                vector<double>          &b,
                vector<double>          &c); 
  //  
  void mTv_mult(vector<vector<double> > &a,
                vector<double>          &b,
                vector<double>          &c); 
  //
  void mv_mult(vector<vector<double> > &a,
                vector<double>          &b,
                vector<double>          &c); 
  //
  void mm_mult(vector<vector<double> > &a,
               vector<vector<double> > &b,
               vector<vector<double> > &c);
  //
  double dprod(vector<double> &a,
               vector<double> &b,
               int            &a_num);
  //
  double dprod(vector<double> &a,
               vector<double> &b);
  //
  void mTm_mult(vector<vector<double> > &a,
                vector<vector<double> > &b,
                vector<vector<double> >       &c);
  //
  void mmT_mult(vector<vector<double> > &a,
                vector<vector<double> > &b,
                vector<vector<double> >       &c);
  //
  void amaT_mult(vector<vector<double> > &a,
                 vector<vector<double> > &b,
                 vector<vector<double> >       &c);
  //
  void aTma_mult(vector<vector<double> > &a,
                 vector<vector<double> > &b,
                 vector<vector<double> >       &c);
  //
  vector<vector<double> > compTransMat(vector<double>& a_x0,
                                       vector<double>& a_x1,
                                       vector<double>& a_x2); 
  //
  vector<double> cross(vector<double>& v1,
                       vector<double>& v2);
  //
  void normalizeVector(vector<double>& u);
  //
  void normalizeVector3(vector<double>& u);
  //
  void wait();
  void wait(string name);
  //
  vector<double> normVec(vector<double>& u);
  //
  void  matrixTranspose(vector<vector<double> >& a,
                        vector<vector<double> >& b);
  //
  double vectorMag(vector<double>& u);
  //
  void dumpMat(vector<vector<double> > &a,string name);
  //
  void dumpMat(vector<vector<int> > &a,string name);
  //
  void dumpVec(vector<double> &a,string name);
  //
  void dumpVec(vector<int> &a,string name);
  //
  bool mat_inverse(vector<vector<double> >&a);
  //
  void dumpVecAsMatMPI(vector<int>& a,int row,string name,int& procID,int&numprocs);
  //
  void dumpVecAsMat(vector<int>& a,int row,string name);
  //
  void dumpVecAsMat(vector<double>& a,int row,string name);
  //
}
#endif
