#define PLEIADES 0
#ifndef KARMA_UTILS_MPI_OPS_H
#define KARMA_UTILS_MPI_OPS_H

#ifdef KARMA_MPI
    #include "mpi.h"
#if(PLEIADES==1)
    #include "mpi_cxx_bindings.h"
#endif
#endif

namespace KARMA {
    
  class MPI_Ops {
    
  public:
    static int Rank     (void);
    static int Size     (void);
    static void mpipause(void);
    
  };

} // end namespace KARMA

#endif
