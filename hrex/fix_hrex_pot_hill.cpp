// clang-format off

#include "fix_hrex_pot_hill.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

FixHREXPotHill::FixHREXPotHill(LAMMPS* lmp, int narg, char** arg)
  : FixHREXPot(lmp, narg, arg)
  , unlambed_pe(0.)
{
  scalar_flag = 1;
  global_freq = 1;
  energy_global_flag = 1;

  // parse user input

  if(narg < 8) {
    error->all(FLERR, "Illegal fix hrex/pot/hill command");
  }

  xflag = yflag = zflag = 1;
  
  if(strcmp(arg[3], "NULL") == 0) xflag = 0;
  else xc = utils::numeric(FLERR, arg[3], false, lmp);
  if(strcmp(arg[4], "NULL") == 0) yflag = 0;
  else yc = utils::numeric(FLERR, arg[4], false, lmp);
  if(strcmp(arg[5], "NULL") == 0) zflag = 0;
  else zc = utils::numeric(FLERR, arg[5], false, lmp);
  h = utils::numeric(FLERR, arg[6], false, lmp);
  sigma = utils::numeric(FLERR, arg[7], false, lmp);
  
  k = 1. / (2. * sigma * sigma);
  r2cut = (4. * sigma) * (4. * sigma);
}

double FixHREXPotHill::root_unlambed_pe()
{
  return unlambed_pe;
}

int FixHREXPotHill::setmask()
{
  return POST_FORCE;
}

void FixHREXPotHill::setup(int vflag)
{
  post_force(vflag);
}

/**
 * Force on atoms due to Gaussian shaped energy penalty
 */
void FixHREXPotHill::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double my_unlambed_pe = 0;
  double lambda = get_lambda();

  for(int i = 0; i < nlocal; i++) {
    if(mask[i] & groupbit) {
      double dx = x[i][0] - xc;
      double dy = x[i][1] - yc;
      double dz = x[i][2] - zc;
      domain->minimum_image(dx, dy, dz);
      if(!xflag) dx = 0;
      if(!yflag) dy = 0;
      if(!zflag) dz = 0;
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 < r2cut) {
	double unlambed_pe_i = h * exp(-k*dr2);
	my_unlambed_pe += unlambed_pe_i;
	double fmag = 2. * k * unlambed_pe_i * lambda;
	f[i][0] += fmag * dx;
	f[i][1] += fmag * dy;
	f[i][2] += fmag * dz;
      }
    }
  }
  MPI_Allreduce(&my_unlambed_pe, &unlambed_pe, 1, MPI_DOUBLE, MPI_SUM, world);
}

double FixHREXPotHill::compute_scalar()
{
  return get_lambda() * unlambed_pe;
}
