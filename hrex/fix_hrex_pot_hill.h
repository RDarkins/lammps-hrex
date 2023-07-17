#ifdef FIX_CLASS
// clang-format off
FixStyle(hrex/pot/hill, FixHREXPotHill);
// clang-format on
#else

#ifndef LMP_FIX_HREX_POT_HILL_H
#define LMP_FIX_HREX_POT_HILL_H

#include "fix_hrex_pot.h"

namespace LAMMPS_NS {

class FixHREXPotHill : public FixHREXPot {
public:
  FixHREXPotHill(class LAMMPS*, int, char**);
  
  int setmask() override;
  void setup(int) override;
  void post_force(int) override;
  double compute_scalar() override;
  double root_unlambed_pe() override;
  
private:
  double unlambed_pe;
  double xc;
  double yc;
  double zc;
  double h;
  double sigma;
  double k;
  double r2cut;
  int xflag;
  int yflag;
  int zflag;
};

} // LAMMPS_NS

#endif
#endif
