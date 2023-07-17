#ifndef LMP_FIX_HREX_H
#define LMP_FIX_HREX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHREX : public Fix {
  
public:
  FixHREX(class LAMMPS*, int, char**);
  
  /**
   * Called after Fix::init() but before the first time step of an hrex command
   */
  virtual void hrex_init() {}

  /**
   * Called at the completion of an hrex command
   */
  virtual void hrex_end() {}

  /**
   * Called immediately after a replica exchange event
   */
  virtual void hrex_swap() {}
  
  int get_ireplica() const;
  int get_nreplicas() const;
  int get_iwalker() const;
  int get_nwalkers() const;
  double get_lambda() const;
  int root_MPI_Allgather_walkers(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, int recvcount, MPI_Datatype recvtype);
  void set_hrex(class HREX*);
  
private:
  class HREX* hrex;
};

} // LAMMPS_NS

#endif
