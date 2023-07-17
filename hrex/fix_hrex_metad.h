#ifndef LMP_FIX_HREX_METAD_H
#define LMP_FIX_HREX_METAD_H

#include "fix_hrex.h"
#include "hrex.h"
#include <vector>

namespace LAMMPS_NS {

typedef std::vector<double> CV;
typedef std::vector<double> Gaussian;

/**
 * Members with 'root_' prefixes will only be invoked by the
 * root (rank 0) MPI process of each partition.
 */
class FixHREXMetad : public FixHREX {
public:
  FixHREXMetad(class LAMMPS* lmp, int narg, char** arg) : FixHREX(lmp, narg, arg) {}
  
  int setmask() final override;
  void end_of_step() final override;

  /**
   * CV for current configuration.
   */
  virtual CV root_collective_variables() = 0;

  /**
   * Bias potential at specified CV.
   */
  virtual double root_bias_potential(CV const&) = 0;

  /**
   * Next Gaussian to be deposited to the bias potential,
   * but do not actually deposit it here.
   */
  virtual Gaussian root_build_gaussian() = 0;

  /**
   * Deposit the specified Gaussian to the bias potential
   */
  virtual void root_deposit_gaussian(Gaussian const&) = 0;

  /**
   * Provide a copy of the entire bias potential packaged
   * into a vector<double>, e.g. a grid representation.
   *
   * This and root_unpack_potential() are used to swap
   * bias potentials between replicas.
   */
  virtual std::vector<double> root_pack_potential() = 0;

  /**
   * Reverse root_pack_potential() so as to replace the
   * current bias potential with that provided.
   */
  virtual void root_unpack_potential(std::vector<double>&) = 0;
};

} // LAMMPS_NS

#endif
