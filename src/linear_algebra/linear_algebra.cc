/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * linear_algebra.cc — Dense and sparse linear algebra routines.
 *   Provides matrix-vector operations used in element assembly and
 *   the built-in parallel conjugate-gradient solver (non-PETSc path).
 */
#include "linear_algebra.h"
#include <iostream>
#include "../fem/options.h"
using namespace std;

namespace KARMA {
  double checksign(double& num)
  {
    double one;
    if (num >= 0.) one = 1.;
    if (num <  0.) one = -1.;
    return one;
  }
  
  vector<double> cross(vector<double>& u,
                       vector<double>& v)
  {   
    //3 dimensional vector cross product
    vector<double> vCross(3,0.);
    vCross[0] = u[1]*v[2] - u[2]*v[1];
    vCross[1] = u[2]*v[0] - u[0]*v[2];
    vCross[2] = u[0]*v[1] - u[1]*v[0];
    return vCross;
  }

  void wait()
  {
#if (DEDICATED_FEM_PROC == 0)
    int mypeno = -1;
    int ierr = MPI_Comm_rank(MPI_COMM_WORLD,&mypeno);
    ierr=MPI_Barrier(MPI_COMM_WORLD);
    if (mypeno==0)
      {
#endif
        cout << "waiting."<<endl;
        int dummy;
        cin >> dummy;
#if (DEDICATED_FEM_PROC == 0)
      }
#endif
    return;
  }
  
  void wait(string name)
  {    
    int mypeno = -1;
    int ierr = MPI_Comm_rank(MPI_COMM_WORLD,&mypeno);
    ierr=MPI_Barrier(MPI_COMM_WORLD);
    if (mypeno==0)
      {
        cout << name << endl;    
        int dummy;
        cin >> dummy;
      }
    return;
  }
  
  void normalizeVector(vector<double>& u)
  {   
    //normalize vector to unit length
    int dim = u.size();
    double mag = 0.;
    for (int i=0;i<dim;++i) { 
      mag += pow(abs(u[i]),2);
    }
    mag=pow(mag,0.5);
    if (mag > 0) { 
      for (int i=0;i<dim;++i) { 
        u[i] /= mag;
      } 
    } else if (mag <= 0) { 
      for (int i=0;i<dim;++i) { 
        u[i] = u[i];
      } 
    }
  }
  
  void normalizeVector3(vector<double>& u)
  {   
    //normalize vectors stored in 1D vector
    int dim = u.size();
    dim /= 3;
    for (int i=0;i<dim;++i)
      {
        double mag = 0.;
        for (int idir=0;idir<3;++idir)
          {
            mag += u[i*3 + idir]*u[i*3 + idir];
          }
        mag = pow(mag,0.5);
        if (mag > 1.E-14)
          {
            for (int idir=0;idir<3;++idir)
              {
                u[i*3 + idir] /= mag;
              } 
          }
        else
          {
            for (int idir=0;idir<3;++idir)
              {
                u[i*3 + idir] = u[i*3 + idir];
              } 
          }
      }

  }
  
  double vectorMag(vector<double>& u)
  {  
    int dim = u.size();
    double mag= 0.;
    for (int i=0;i<dim;++i)
      {
	mag += pow(abs(u[i]),2);
      }
    mag=pow(mag,0.5);
    return mag;
  }

  vector<double> normVec(vector<double>& u)
  {   

    double sum=0.;
    for (int v=0;v<u.size();v++) {
      sum+=u[v]*u[v];
    }
    double sumsqrt=pow(sum,0.5);
    vector<double> normVec(u.size(),0.);
    for (int v=0;v<u.size();v++) normVec[v]=u[v]/sumsqrt;
    return normVec;
  }
  
  void mTv_mult(vector<vector<double> > &a,
                vector<double> &b,
                vector<double> &c) 
  {//c = a*b (a=mat, b=vec, c=vec)

    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int mc=c.size();

    for (int m=0;m<ma;++m) 
      {
        c[m]=0.;          
        for (int n=0;n<na;++n) 
          {
            c[m]+=a[m][n]*b[n];
          }
      }
  }
  
  void mv_mult(vector<vector<double> > &a,
               vector<double> &b,
               vector<double> &c)
  {

    fill(c.begin(),c.end(),0.);
    const int nr=a.size();
    const int nc=a[0].size();

    for (int r=0;r<nr;++r)
      {
        for (int cc=0;cc<nc;++cc)
          {
            c[r] += a[r][cc]*b[cc];
          }
      }

  }
  
  void mm_mult(std::vector<std::vector<double> > &a,
               std::vector<std::vector<double> > &b,
               std::vector<std::vector<double> > &c) 
  {//c=a*b
    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int nb=b[0].size();
    const int mc=c.size();
    const int nc=c[0].size();
    for (int m=0;m<ma;++m) 
      {
        for (int n=0;n<nb;++n) 
          {
            c[m][n]=0.;
            //#pragma unroll
            for (int k=0;k<na;++k) 
              {
                c[m][n]+=a[m][k]*b[k][n];
              }
          }
      }
  }
  
  double dprod(vector<double> &a,
               vector<double> &b,
               int&          a_num)
  {
    double dotp=0.;          
    for (int m=0;m<a_num;++m) 
      {
        dotp += a[m]*b[m];
      }
    return dotp;  
  }

  double dprod(vector<double> &a,
               vector<double> &b)
  {
    double dotp=0.;          
    for (int m=0;m<a.size();++m) 
      {
        dotp += a[m]*b[m];
      }
    return dotp;  
  }
  
  void mTm_mult(vector<vector<double> > &a,
                vector<vector<double> > &b,
                vector<vector<double> > &c) 
  {//c=transpose(a)*b
    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int nb=b[0].size();
    const int mc=c.size();
    const int nc=c[0].size();
    
    for (int m=0;m<na;++m) 
      {
        for (int n=0;n<nb;++n) 
          {
            c[m][n]=0.;
            for (int k=0;k<ma;++k) 
              {
                c[m][n]+=a[k][m]*b[k][n];
              }
          }
      }
  }
  
  void mmT_mult(vector<vector<double> > &a,
                vector<vector<double> > &b,
                vector<vector<double> > &c)
  {//c=a*transpose(b)
    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int nb=b[0].size();
    const int mc=c.size();
    const int nc=c[0].size();
    
    for (int m=0;m<ma;++m) 
      {
        for (int n=0;n<mb;++n) 
          {
            c[m][n]=0.;
            for (int k=0;k<na;++k) 
              {
                c[m][n]+=a[m][k]*b[n][k];
              }
          }
      }
  }
  
  void amaT_mult(vector<vector<double> > &a,
                 vector<vector<double> > &b,
                 vector<vector<double> > &c)
  {//c=a*b*transpose(a)
    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int nb=b[0].size();
    const int mc=c.size();
    const int nc=c[0].size();
    
    vector<double>       dummy(nb,0.0);
    vector<vector<double> > ab(ma,dummy);
    mm_mult(a,b,ab);
    
    mmT_mult(ab,a,c);
    
  }
  
  void aTma_mult(vector<vector<double> > &a,
                 vector<vector<double> > &b,
                 vector<vector<double> > &c)
  {//c=transpose(a)*b*a
    const int ma=a.size();
    const int na=a[0].size();
    const int mb=b.size();
    const int nb=b[0].size();
    const int mc=c.size();
    const int nc=c[0].size();

    vector<double>       dummy(na,0.0);
    vector<vector<double> > ba(mb,dummy);
    mm_mult(b,a,ba);
    mTm_mult(ba,a,c);  
  }
  
  vector<vector<double> > compTransMat(vector<double>& a_x0,
				       vector<double>& a_x1,
				       vector<double>& a_x2)
  {
    vector<double> v1(a_x0.size(),0.);
    vector<double> v2(a_x0.size(),0.);

    for (int i=0;i<a_x1.size();++i) {
      v1[i]=a_x1[i]-a_x0[i];
      v2[i]=a_x2[i]-a_x0[i];
    }
    vector<double> v3(a_x0.size(),0.);
    v3=cross(v1,v2);
    v2=cross(v3,v1);
    vector<double> e1(a_x0.size(),0.);
    vector<double> e2(a_x0.size(),0.);
    vector<double> e3(a_x0.size(),0.);
    e1 = normVec(v1);
    e2 = normVec(v2);
    e3 = normVec(v3);

    vector<double> dummy(3,0.);
    vector<vector<double> > T(3,dummy);
    for (int v=0;v<3;++v) {
      T[0][v]=e1[v]; 
      T[1][v]=e2[v];
      T[2][v]=e3[v];
    }
    return T;
  }

  void matrixTranspose(vector<vector<double> >& a,
                       vector<vector<double> >& b) 
  {
    for (int i=0;i<a.size();++i) {
      for (int j=0;j<a[0].size();++j) {
        b[j][i] = a[i][j];
      }
    }
  }

  void dumpMat(vector<vector<double> >& a,string name)
  {
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"
    #define GREEN "\033[1m\033[32m"
    cout << BLUE << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;
    for (int j=0;j<a.size();++j) {
      for (int i=0;i<a[0].size();++i) {
        if (a[j][i] == 0.) {
          cout << a[j][i] << " ";
        } else {
          cout << GREEN << a[j][i] <<  " " << RESET;
        }
      }
      cout <<endl;  
    }
  }

  void dumpMat(vector<vector<int> >& a,string name)
  {
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"
    #define GREEN "\033[1m\033[32m"
    cout << BLUE << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;
    for (int j=0;j<a.size();++j) {
      for (int i=0;i<a[0].size();++i) {
        if (a[j][i] == 0.) {
          cout << a[j][i] << " ";
        } else {
          cout << GREEN << a[j][i] <<  " " << RESET;
        }
      }
      cout <<endl;  
    }
  }
  
  void dumpVec(vector<double>& a,string name)
  {
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"
    cout << BLUE << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;
    for (int j=0;j<a.size();++j) {
      if (a[j] == 0.) {
        cout << a[j] << " ";
      } else {
        cout << GREEN << a[j] <<  " " << RESET;
      }
    }
    cout << endl;
  }

  void dumpVec(vector<int>& a,string name)
  {
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"
    cout << BLUE << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;
    for (int j=0;j<a.size();++j) {
      if (a[j] == 0.) {
        cout << a[j] << " ";
      } else {
        cout << GREEN << a[j] <<  " " << RESET;
      }
    }
    cout << endl;
  }
  
  void dumpVecAsMat(vector<double>& a,int row,string name,int& procID,int&numprocs)
  {
    //colors
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"

    int cutoff = a.size()/row;
    int ierr;

    for (int np=0;np<numprocs;++np)
      {
        ierr = MPI_Barrier(MPI_COMM_WORLD);            
        if (np == procID)
          {

            cout << BLUE << procID << " | " << name  << ": "<< endl << RESET;
            cout.precision(2);
            cout << scientific;
            cout << std::showpos;

            for (int j=0;j<cutoff;++j)
              {
                for (int i=0;i<row;++i)
                  {
                    int index = i + j*row;
                    if (abs(a[index]) == 0)
                      {
                        cout << a[index] << " ";
                      }
                    else
                      {
                        cout << GREEN << a[index] <<  " " << RESET;
                      }
                  }
                cout << endl;
              }
          }
      }
  }

  void dumpVecAsMatMPI(vector<double>& a,int row,string name,int& procID,int&numprocs)
  {
    //colors
    cout << endl;
    #define BLUE "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"

    int cutoff = a.size()/row;
    int ierr;

    for (int np=0;np<numprocs;++np)
      {
        ierr = MPI_Barrier(MPI_COMM_WORLD);            
        if (np == procID)
          {

            cout << BLUE << procID << " | " << name  << ": "<< endl << RESET;
            cout.precision(2);
            cout << scientific;
            cout << std::showpos;

            for (int j=0;j<cutoff;++j)
              {
                for (int i=0;i<row;++i)
                  {
                    int index = i + j*row;
                    if (abs(a[index]) == 0)
                      {
                        cout << a[index] << " ";
                      }
                    else
                      {
                        cout << GREEN << a[index] <<  " " << RESET;
                      }
                  }
                cout << endl;
              }
          }
      }
  }

  void dumpVecAsMat(vector<int>& a,int row,string name)
  {
    //colors
    cout << endl;
    #define BLUE  "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"

    int cutoff = a.size()/row;
    int ierr;

    cout << BLUE << " | " << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;

    for (int j=0;j<cutoff;++j)
      {
        for (int i=0;i<row;++i)
          {
            int index = i + j*row;
            if (abs(a[index]) < 1.E-20)
              {
                cout << a[index] << " ";
              }
            else
              {
                cout << GREEN << a[index] <<  " " << RESET;
              }
          }
        cout << endl;
      }

  }

  void dumpVecAsMat(vector<double>& a,int row,string name)
  {
    //colors
    cout << endl;
    #define BLUE  "\033[1m\033[31m" 
    #define RESET "\033[0m"  
    #define GREEN "\033[1m\033[32m"

    int cutoff = a.size()/row;
    int ierr;

    cout << BLUE << " | " << name  << ": "<< endl << RESET;
    cout.precision(2);
    cout << scientific;
    cout << std::showpos;

    for (int j=0;j<cutoff;++j)
      {
        for (int i=0;i<row;++i)
          {
            int index = i + j*row;
            if (abs(a[index]) < 1.E-20)
              {
                cout << a[index] << " ";
              }
            else
              {
                cout << GREEN << a[index] <<  " " << RESET;
              }
          }
        cout << endl;
      }

  }

  bool mat_inverse(vector<vector<double> > &a) 
  {    
    int n=a.size();
    double pvt;
    vector<double> tmpVec(n,0.);
    double tmp;
    int i,j,k;
    double EPS=1.e-12;

    vector<vector<double> > b(n,vector<double>(a[0].size(),0.));
    for (int d=0;d<n;++d) b[d][d]=1.;

    for (i=0;i<n;i++) 
      {
        //get pivot value
        pvt = a[i][i];

        if (abs(pvt) < EPS) 
          {
            for (j=i+1;j<n;++j) 
              {
                if(abs(pvt = a[j][i]) >= EPS) break;
              }
            if (abs(pvt) < EPS) 
              {
                cout <<"MATRIX ill conditioned" << "\t" << pvt << endl;
                return false;     // nowhere to run!
              }
            a[i].swap(a[j]);
            tmpVec=b[j];
            b[j]=b[i];
            b[i]=tmpVec;
          }
        //Gaussian elimination of column
        for (k=i+1;k<n;++k) 
          {
            tmp = a[k][i]/pvt;

            for (j=i+1;j<n;++j) 
              {
                a[k][j] -= tmp*a[i][j];
              }
            for (j=0;j<n;++j) 
              {
                b[k][j] -= tmp*b[i][j];
              }
          }
      } 
    //back substitution
    for (i=n-1;i>=0;i--) 
      {
        for (k=0;k<n;++k) 
          {
            for (j=n-1;j>i;j--) 
              {
                b[i][k] -= a[i][j]*b[j][k];
              }
            b[i][k] /= a[i][i];
          }
      }    

    for (int r=0;r<n;++r)
      {
        for (int c=0;c<a[0].size();++c)
          {
            a[r][c] = b[r][c];
          }
      }

    return true;
  }
  
} // end namespace KARMA

