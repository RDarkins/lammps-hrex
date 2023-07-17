#ifndef LMP_FIX_HREX_DUMP_H
#define LMP_FIX_HREX_DUMP_H

#include "fix_hrex.h"

namespace LAMMPS_NS {

class FixHREXDump : public FixHREX {
public:
  FixHREXDump(class LAMMPS* lmp, int narg, char** arg) : FixHREX(lmp, narg, arg) {}
  
  /**
   * Called at the start of an hrex command before both Fix::init and FixHREX::hrex_init.
   * ireplica will be replaced with the returned value. Useful for restarts.
   */
  virtual int restore_ireplica() {return get_ireplica();}

  virtual void hrex_init() override {}
  virtual void hrex_end() override {}
  virtual void hrex_swap() override {}
  virtual int setmask() override {return 0;}
};

} // LAMMPS_NS

#endif
