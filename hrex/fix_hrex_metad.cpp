#include "fix_hrex_metad.h"

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

int FixHREXMetad::setmask() 
{
  return POST_FORCE | END_OF_STEP;
}

/**
 * Multiple walkers are implemented here so deriving
 * classes do not need to deal with them.
 */
void FixHREXMetad::end_of_step()
{
  if(comm->me) return;

  int nwalkers = get_nwalkers();
  Gaussian g = root_build_gaussian();
  if(nwalkers == 1) {
    root_deposit_gaussian(g);
  } else {
    // gather all gaussians from walkers
    
    int gsize = g.size();
    std::vector<double> all_gaussians(gsize * nwalkers);
    root_MPI_Allgather_walkers(g.data(), gsize, MPI_DOUBLE,
			       all_gaussians.data(), gsize, MPI_DOUBLE);
    
    // deposit all gaussians

    for(int i=0; i<nwalkers; i++) {
      std::vector<double> gi(all_gaussians.begin() + i * gsize,
			     all_gaussians.begin() + (i + 1) * gsize);
      root_deposit_gaussian(gi);
    }
  }
}
