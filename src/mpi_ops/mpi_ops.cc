#include "mpi_ops.h"

using namespace std;

#ifdef KARMA_MPI
using namespace MPI;
#endif

namespace KARMA {
   
  int MPI_Ops::Rank(void) {
#ifdef KARMA_MPI
    int output; 
    MPI_Comm_rank(MPI_COMM_WORLD,&output);
    return output;
#else 
    return 0
#endif
      }
  
  int MPI_Ops::Size(void) {
    int output; 
    MPI_Comm_size(MPI_COMM_WORLD,&output);
    return output;

  }
} // end namespace KARMA
