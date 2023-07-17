#ifdef FIX_CLASS
// clang-format off
FixStyle(hrex/dump/swap, FixHREXDumpSwap);
// clang-format on
#else

#ifndef LMP_FIX_HREX_DUMP_SWAP_H
#define LMP_FIX_HREX_DUMP_SWAP_H

#include "fix_hrex_dump.h"

namespace LAMMPS_NS {

class FixHREXDumpSwap : public FixHREXDump {
public:
  FixHREXDumpSwap(class LAMMPS*, int, char**);
  ~FixHREXDumpSwap() override;
  
  int restore_ireplica() override;
  void hrex_init() override;
  void hrex_end() override;
  void hrex_swap() override;
  
private:
  FILE* fp;
  char* filename;
};

} // LAMMPS_NS

#endif
#endif
