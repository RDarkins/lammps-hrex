#ifndef LMP_FIX_HREX_POT_H
#define LMP_FIX_HREX_POT_H

#include "fix_hrex.h"

namespace LAMMPS_NS {

/**
 * Members with 'root_' prefixes will only be invoked by the
 * root (rank 0) MPI process of each partition.
 */
class FixHREXPot : public FixHREX {
public:
  FixHREXPot(class LAMMPS* lmp, int narg, char** arg) : FixHREX(lmp, narg, arg) {}
  
  /**
   * Return the potential energy modification dU for the current
   * configuration. It must not be multiplied by lambda_i.
   */
  virtual double root_unlambed_pe() = 0;
};

} // LAMMPS_NS

#endif
