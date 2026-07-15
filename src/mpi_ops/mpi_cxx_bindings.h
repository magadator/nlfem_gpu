#ifndef KARMA_MPI_CXX_BINDINGS_H
#define KARMA_MPI_CXX_BINDINGS_H 
#include "mpi.h"
namespace MPI
{
  const static decltype(MPI_INT) INT=MPI_INT;
  const static decltype(MPI_INT) INTEGER=MPI_INT;  
  const static decltype(MPI_INT) DOUBLE=MPI_DOUBLE;
  const static decltype(MPI_INT) BOOL=MPI_LOGICAL;
}
#endif
