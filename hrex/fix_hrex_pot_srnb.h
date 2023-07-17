#ifdef FIX_CLASS
// clang-format off
FixStyle(hrex/pot/srnb, FixHREXPotSRNB);
// clang-format on
#else

#ifndef LMP_FIX_HREX_POT_SRNB_H
#define LMP_FIX_HREX_POT_SRNB_H

#include "fix_hrex_pot.h"

namespace LAMMPS_NS {

class FixHREXPotSRNB : public FixHREXPot {
public:
  FixHREXPotSRNB(class LAMMPS*, int, char**);
  
  int setmask() override;
  void init() override;
  void init_list(int, NeighList*) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  double compute_scalar() override;
  double root_unlambed_pe() override;
  
private:
  int jgroup, jgroupbit;
  double scale_factor;
  double unlambed_pe;
  class NeighList* list;
};

} // LAMMPS_NS

#endif
#endif
